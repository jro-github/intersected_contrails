#!/bin/bash

# === 用户参数配置 ===
GITHUB_REPO_URL="https://github.com/MingchuanLuo/intersected_contrails.git"  # 修改为你的仓库地址
BRANCH_NAME="main"                     # main 或 master
NEW_TAG="v2.0.0"                       # 可选的版本标签（如不需要请留空）

# === 检查依赖 ===
echo "🔍 检查是否安装 Git..."
if ! command -v git &> /dev/null; then
  echo "❌ Git 未安装，请先安装 Git。"
  exit 1
fi

# === 开始操作 ===
echo "🧹 清理旧的 Git 信息（如果有）..."
rm -rf .git

echo "🛠️ 初始化新的 Git 仓库..."
git init
git remote add origin "$GITHUB_REPO_URL"

echo "📦 添加所有新文件..."
git add .

echo "✅ 提交新代码..."
git commit -m "Replace project with new version"

echo "🚀 强制推送到远程仓库分支：$BRANCH_NAME ..."
git branch -M "$BRANCH_NAME"
git push origin "$BRANCH_NAME" --force

# === 打标签 ===
if [[ -n "$NEW_TAG" ]]; then
    echo "🏷️ 打版本标签：$NEW_TAG"
    git tag "$NEW_TAG"
    git push origin "$NEW_TAG"
fi

echo "🎉 更新完成！请在 GitHub 上检查你的仓库。"

